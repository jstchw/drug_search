import { useUrlParams } from "../../hooks/useUrlParams";
import { useNavigate } from "react-router-dom";
import SearchBox from "../../components/SearchBox/SearchBox";
import Header from "../../components/Header/Header";
import { Col, Row } from "react-bootstrap";
import useDrugInfo from "../../hooks/useDrugInfo";
import React from "react";
import DrugPropertyBox from "../../components/DrugPropertyBox/DrugPropertyBox";
import { SearchHistoryContext } from "../../contexts/SearchHistoryContext";
import ChartSection from "../../components/ChartSection/ChartSection";
import {DrugProperties} from "../../types";
import {Badge} from "react-bootstrap";

const Result = () => {
    const { params, error } = useUrlParams()
    const capitalizedTerms = params.terms.map(term => term.charAt(0).toUpperCase() + term.slice(1))
    const navigate = useNavigate()

    // Retrieve drug info from the API (DrugSearch server)
    const drugInfo = useDrugInfo(params)

    const groupedByBrandName = React.useMemo(() => {
        if (params.searchBy === 'brand_name') {
            return drugInfo.reduce<Record<string, DrugProperties[]>>((acc, drug) => {
                if (drug.product) {
                    (acc[drug.product] = acc[drug.product] || []).push(drug)
                }
                return acc
            }, {})
        } else {
            return drugInfo
        }
    }, [drugInfo, params.searchBy])

    console.log(groupedByBrandName)

    // Effect to redirect to the error page if there is an error
    React.useEffect(() => {
        if (error) {
            navigate('/error')
        }
    }, [error, navigate])

    // Effect to update the document title
    React.useEffect(() => {
        if (params.terms) {
            document.title = capitalizedTerms.join(' & ') + ' - Drug Search'
        }
    }, [capitalizedTerms, params.terms])

    // Updating search history (setting only)
    const searchHistoryContext = React.useContext(SearchHistoryContext)
    const { updateSearchHistory } = searchHistoryContext || {};

    React.useEffect(() => {
        if(drugInfo.length > 0 && updateSearchHistory) {
            updateSearchHistory(params)
        }
    }, [drugInfo, params, updateSearchHistory])

    return (
        <>
            <div className="d-flex flex-column justify-content-center align-items-center">
                <Row className={'mb-4 mt-5'}>
                    <Header/>
                </Row>
                <Row className={'mb-4 text-center'}>
                    <SearchBox passedInput={capitalizedTerms}/>
                </Row>
            </div>

            <Row className={'justify-content-center mx-auto'}>
                {params.searchBy === 'brand_name' ? (
                    Object.entries(groupedByBrandName).map(([brandName, drugs], index) => (
                        <React.Fragment key={index}>
                            <Col xs={12} className={'mb-1 d-flex justify-content-center'}>
                                <h1>{brandName}</h1>
                            </Col>
                            <Row className={'d-flex justify-content-center mb-2'}>
                                <Badge className={'fst-italic'} style={{width: 'fit-content'}}>Active substances</Badge>
                            </Row>
                            {(drugs as DrugProperties[]).map((drug, index) => (
                                <Col xs={drugs.length === 1 ? 6 : 4} key={index} className="mb-4">
                                    <DrugPropertyBox drug={drug} isSingle={drugs.length === 1}/>
                                </Col>
                            ))}
                        </React.Fragment>
                    ))
                ) : (
                    drugInfo.map((drug, index) => (
                        <Col xs={drugInfo.length === 1 ? 6 : 4} key={index} className="mb-4">
                            <DrugPropertyBox drug={drug} isSingle={drugInfo.length === 1}/>
                        </Col>
                    ))
                )}
            </Row>

            <Row className={'justify-content-center mx-auto'}>
                <Col xs={6} className="mb-4">
                    <ChartSection />
                </Col>
            </Row>
        </>
    );
}

export default Result;