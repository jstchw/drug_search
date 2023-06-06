import React from "react";
import {Col, Container, Row, Badge, Placeholder, Accordion} from "react-bootstrap";
import Histogram from "../Histogram/Histogram";
import './SearchResultObject.css'
import DrugDescription from "../DrugDescription/DrugDescription";
import AccordionBody from "react-bootstrap/AccordionBody";


const SearchResultObject = ({ searchResults, searchTerm }) => {
    const { termCountDict, totalCount } = searchResults.result
    const [retrievedTermArr, setRetrievedTermArr] = React.useState(null)

    const processDrugGroups = (groups) => {
        return groups.map((group, index) => {
            let variant
            switch (group) {
                case 'approved':
                    variant = 'success'
                    break
                case 'investigational':
                    variant = 'warning'
                    break
                case 'illicit':
                    variant = 'danger'
                    break
                case 'vet_approved':
                    variant = 'info'
                    break
                default:
                    variant = 'secondary'
            }
            group = group.substring(0, 3).toUpperCase()
            return <Badge className={'mx-1'} key={index} bg={variant}>{group}</Badge>
        })

    }

    return (
        <div>
            <Container className="my-4">
                {retrievedTermArr ? (
                    <>
                        {retrievedTermArr && <h1>{retrievedTermArr[0]}</h1>}
                        {retrievedTermArr && <h3>{processDrugGroups(retrievedTermArr[1])}</h3>}
                    </>
                ) : (
                    <>
                        <Placeholder as={'h1'} animation="glow">
                            <Placeholder xs={6} />
                        </Placeholder>
                        <Placeholder as={'h3'} animation="glow">
                            <Placeholder xs={3} />
                        </Placeholder>
                    </>
                )}
            </Container>
            <Container className={'main-drug-container'}>
                <Row>
                    <Col className={'alert alert-secondary mx-2'}>
                        <span>The provided information is indicative and shouldn't be used for inference.</span>
                    </Col>
                </Row>
                <Row>
                    <Col>
                        <DrugDescription drugName={searchTerm} onRetrieved={setRetrievedTermArr}/>
                    </Col>
                    <Col>
                        <Accordion>
                            <Accordion.Item eventKey={'0'}>
                                <Accordion.Header><Badge>Half-life</Badge></Accordion.Header>
                                <Accordion.Body>
                                    {retrievedTermArr ? (
                                        retrievedTermArr[4]
                                    ) : (
                                        <Placeholder animation="glow">
                                            <Placeholder xs={7} /> <Placeholder xs={4} /> <Placeholder xs={4} />{' '}
                                            <Placeholder xs={6} /> <Placeholder xs={8} />
                                        </Placeholder>
                                    )}
                                </Accordion.Body>
                            </Accordion.Item>
                        </Accordion>
                    </Col>
                </Row>
                <Row>
                    <Col className={'mt-3'}>
                        {searchResults && <Histogram termCountDict={termCountDict} totalCount={totalCount} type={'all_groups'} />}
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default SearchResultObject;