import React, {useEffect} from 'react';
import {Modal, Form, InputGroup, ToggleButton, Button} from 'react-bootstrap';
import countries from '../../assets/countries.json';
import {ThemeContext} from "../../contexts/ThemeContext";
import {Moon, Sun} from "react-bootstrap-icons";
import Cookies from 'js-cookie'


// Search types that can be selected in the popover
export const searchTypes = [
    {
        value: 'patient.drug.openfda.generic_name',
        index: 0,
        label: 'Generic Name',
        type: 'searchBy'
    },
    {
        value: 'patient.drug.openfda.brand_name',
        index: 1,
        label: 'Brand Name',
        type: 'searchBy'
    },
    {
        value: 'patient.reaction.reactionmeddrapt',
        index: 2,
        label: 'Side Effect',
        type: 'searchBy'
    },
]

export const searchSex = [
    {
        value: 'patient.patientsex:1',
        index: 0,
        label: 'Male',
        type: 'sex'
    },
    {
        value: 'patient.patientsex:2',
        index: 1,
        label: 'Female',
        type: 'sex'
    }
]

export const searchAgeRange = [
    {
        value: '',
        index: 0,
        type: 'age',
    },
    {
        value: '',
        index: 1,
        type: 'age',
    }
]

export const searchCountry = Object.entries(countries).map(([value, label], index) => ({
    value,
    label,
    index,
    type: 'country'
}))

const OptionModal = (props) => {
    const isSexDisabledCookie = Cookies.get('isSexDisabled')
    const [isSexDisabled, setIsSexDisabled] = React.useState(
        isSexDisabledCookie === 'false' ? false : isSexDisabledCookie === 'true' ? true : true
    )
    const [lastSelectedSex, setLastSelectedSex] = React.useState(
        Cookies.get('lastSelectedSex') || searchSex[0].value
    )

    const isAgeDisabledCookie = Cookies.get('isAgeDisabled')
    const [isAgeDisabled, setIsAgeDisabled] = React.useState(
        isAgeDisabledCookie === 'false' ? false : isAgeDisabledCookie === 'true' ? true : true
    )

    const ageCookieStart = Cookies.get('lastSelectedAgeStart')
    const ageCookieEnd = Cookies.get('lastSelectedAgeEnd')
    const [ageRange, setAgeRange] = React.useState(
        [ageCookieStart, ageCookieEnd] || [searchAgeRange[0].value, searchAgeRange[1].value]
    )

    const [lastSelectedAge, setLastSelectedAge] = React.useState(
        [ageCookieStart, ageCookieEnd] || [searchAgeRange[0].value, searchAgeRange[1].value]
    )

    const isCountryDisabledCookie = Cookies.get('isCountryDisabled')
    const [isCountryDisabled, setIsCountryDisabled] = React.useState(
        isCountryDisabledCookie === 'false' ? false : isCountryDisabledCookie === 'true' ? true : true
    )
    const defaultCountry = searchCountry.find(country => country.value === 'US')

    const [lastSelectedCountry, setLastSelectedCountry] = React.useState(
        Cookies.get('lastSelectedCountry') || defaultCountry.value
    )

    const { theme, toggleTheme } = React.useContext(ThemeContext)

    useEffect(() => {
        Cookies.set('isSexDisabled', isSexDisabled, { expires: 365 }, { secure: true })
        Cookies.set('isAgeDisabled', isAgeDisabled, { expires: 365 }, { secure: true })
        Cookies.set('isCountryDisabled', isCountryDisabled, { expires: 365 }, { secure: true })
    }, [isSexDisabled, isAgeDisabled, isCountryDisabled])

    // Clear the age range when the age is disabled
    useEffect(() => {
        if(isAgeDisabled) {
            props.setSearchOptions({
                ...props.searchOptions,
                age: undefined
            })
        } else {
            props.setSearchOptions({
                ...props.searchOptions,
                age: lastSelectedAge
            })
        }
    }, [isAgeDisabled]) // eslint-disable-line react-hooks/exhaustive-deps


    useEffect(() => {
        let newSearchOptions = { ...props.searchOptions };

        if(isCountryDisabled) {
            newSearchOptions.country = undefined;
        } else {
            newSearchOptions.country = lastSelectedCountry;
        }

        if(isSexDisabled) {
            newSearchOptions.sex = undefined;
        } else {
            newSearchOptions.sex = lastSelectedSex;
        }

        props.setSearchOptions(newSearchOptions);
    }, [isCountryDisabled, isSexDisabled]);

    const sanitizeAge = (age) => {
        const lastChar = age.slice(-1); // Get the last character
        if(!/^[0-9]*$/.test(lastChar)) {
            age = age.slice(0, -1); // Remove the last character if it's not a digit
        }
        return age;
    }
    const handleOptionsChange = (index, type, value) => {
        switch (type) {
            case 'searchBy':
                props.setSearchOptions({
                    ...props.searchOptions,
                    [type]: searchTypes[index].value
                })
                Cookies.set('searchBy', searchTypes[index].value, { expires: 365 }, { secure: true })

                props.setSelectedSearchOptionIndex({
                    ...props.selectedSearchOptionIndex,
                    [type]: index
                })


                break
            case 'sex':
                setLastSelectedSex(searchSex[index].value)
                Cookies.set('lastSelectedSex', searchSex[index].value, { expires: 365 }, { secure: true })
                props.setSearchOptions({
                    ...props.searchOptions,
                    [type]: searchSex[index].value
                })

                props.setSelectedSearchOptionIndex({
                    ...props.selectedSearchOptionIndex,
                    [type]: index
                })
                break
            case 'age':
                const newAgeRange = [...ageRange]
                newAgeRange[index] = value
                setAgeRange(newAgeRange)
                Cookies.set('lastSelectedAgeStart', newAgeRange[0], { expires: 365 }, { secure: true })
                Cookies.set('lastSelectedAgeEnd', newAgeRange[1], { expires: 365 }, { secure: true })
                if(!isAgeDisabled) {
                    setLastSelectedAge(newAgeRange)
                    props.setSearchOptions({
                        ...props.searchOptions,
                        [type]: newAgeRange
                    })
                }
                break
            case 'country':
                setLastSelectedCountry(searchCountry[index].value)
                Cookies.set('lastSelectedCountry', searchCountry[index].value, { expires: 365 }, { secure: true })
                props.setSearchOptions({
                    ...props.searchOptions,
                    [type]: searchCountry[index].value
                })

                props.setSelectedSearchOptionIndex({
                    ...props.selectedSearchOptionIndex,
                    [type]: index
                })
                break
            default:
                break
        }
    }

    return (
        <>
            <Modal centered show={props.show} onHide={props.handleClose}>
                <Modal.Header closeButton>
                    <Button size={'md'} onClick={toggleTheme} className={'me-2'}>
                        {theme === 'light' ? <Moon /> : <Sun />}
                    </Button>
                    <Modal.Title>Options</Modal.Title>
                </Modal.Header>
                <Modal.Body>
                    <Form>
                        <Form.Group className="mb-3 d-flex align-items-center">
                            <div className={'d-flex align-items-center'}>
                                <ToggleButton
                                    type="checkbox"
                                    variant="outline-primary"
                                    checked={true}
                                    value="1"
                                    disabled={true}
                                >
                                    Search by
                                </ToggleButton>
                            </div>
                            <InputGroup className={'flex-grow-1 mx-3'} style={{width: 'auto'}}>
                                <Form.Select
                                    onChange={(e) => handleOptionsChange(e.target.value, 'searchBy')}
                                    value={props.selectedSearchOptionIndex.searchBy}
                                    style={{width: 'auto'}}
                                >
                                    {searchTypes.map((searchType, index) => (
                                        <option
                                            key={index}
                                            value={searchType.index}
                                        >
                                            {searchType.label}
                                        </option>
                                    ))}
                                </Form.Select>
                            </InputGroup>
                        </Form.Group>

                        <Form.Group className="mb-3 d-flex align-items-center">
                            <div className={'d-flex align-items-center'}>
                                <ToggleButton
                                    type="checkbox"
                                    variant="outline-primary"
                                    checked={!isSexDisabled}
                                    value="1"
                                    onClick={() => setIsSexDisabled(!isSexDisabled)}
                                >
                                    Sex
                                </ToggleButton>
                            </div>
                            <InputGroup className={'mx-3 flex-grow-1'}>
                                <Form.Select
                                    onChange={e => {handleOptionsChange(e.target.value, 'sex')}}
                                    value={props.selectedSearchOptionIndex.sex}
                                    disabled={isSexDisabled}
                                >
                                    {searchSex.map((sex, index) => (
                                            <option
                                                key={index}
                                                value={sex.index}
                                            >
                                                {sex.label}
                                            </option>
                                    ))}
                                </Form.Select>
                            </InputGroup>
                        </Form.Group>

                        <Form.Group className="mb-3 d-flex align-items-center">
                            <div className={'d-flex align-items-center'}>
                                <ToggleButton
                                    type="checkbox"
                                    variant="outline-primary"
                                    checked={!isAgeDisabled}
                                    value="1"
                                    onClick={() => setIsAgeDisabled(!isAgeDisabled)}
                                >
                                Age
                                </ToggleButton>
                            </div>

                            <InputGroup className={'mx-3 flex-grow-1'}>
                                <Form.Control
                                    type="text"
                                    placeholder="Min"
                                    value={ageRange[0]}
                                    onChange={(e) => {
                                        e.target.value = sanitizeAge(e.target.value)
                                        handleOptionsChange(0, 'age', e.target.value);
                                    }}
                                    disabled={isAgeDisabled}
                                    min={0}
                                    onCopy={(e) => (e.preventDefault())}
                                    onPaste={(e) => (e.preventDefault())}
                                />
                                <InputGroup.Text>-</InputGroup.Text>
                                <Form.Control
                                    type="text"
                                    placeholder="Max"
                                    value={ageRange[1]}
                                    onChange={(e) => {
                                        e.target.value = sanitizeAge(e.target.value)
                                        handleOptionsChange(1, 'age', e.target.value);
                                    }}
                                    disabled={isAgeDisabled}
                                    min={0}
                                    onCopy={(e) => (e.preventDefault())}
                                    onPaste={(e) => (e.preventDefault())}
                                />
                            </InputGroup>
                        </Form.Group>
                        <Form.Group>
                            <div className={'d-flex align-items-center'}>
                                <ToggleButton
                                    type="checkbox"
                                    variant="outline-primary"
                                    value="1"
                                    onClick={() => setIsCountryDisabled(!isCountryDisabled)}
                                    checked={!isCountryDisabled}
                                >
                                    Country
                                </ToggleButton>
                                <InputGroup className={'mx-3 flex-grow-1'}>
                                    <Form.Select
                                        disabled={isCountryDisabled}
                                        onChange={(e) => handleOptionsChange(e.target.value, 'country')}
                                        value={props.selectedSearchOptionIndex.country}
                                    >
                                        {searchCountry.map((country, index) => (
                                            <option key={index} value={index}>
                                                {country.label}
                                            </option>
                                        ))}
                                    </Form.Select>
                                </InputGroup>
                            </div>
                        </Form.Group>
                    </Form>
                </Modal.Body>
            </Modal>
        </>
    )
}

export default OptionModal;